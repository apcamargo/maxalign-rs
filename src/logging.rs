//! Logging initialization for the CLI.

use std::fmt;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use nu_ansi_term::{Color, Style};
use thiserror::Error;
use time::{OffsetDateTime, UtcOffset};
use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::fmt::{
    FmtContext,
    format::{FormatEvent, FormatFields, Writer},
};
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::prelude::*;
use tracing_subscriber::registry::LookupSpan;

/// Errors raised while installing global logging subscribers.
#[derive(Debug, Error)]
pub enum InitError {
    #[error("failed to create log file '{path}': {source}")]
    CreateLogFile {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to install tracing subscriber: {0}")]
    SetGlobalDefault(#[from] tracing::subscriber::SetGlobalDefaultError),
}

/// Installs the tracing subscriber for the CLI process.
///
/// This CLI intentionally captures only native `tracing` events. If a future
/// dependency emits important `log` records, that bridge should be added back
/// explicitly rather than carried preemptively.
pub fn init(verbose: u8, log_path: Option<&Path>) -> Result<(), InitError> {
    let clock = TimestampClock::detect();
    let stderr_ansi_enabled = std::io::stderr().is_terminal();
    if let Some(path) = log_path {
        let file = std::fs::File::create(path).map_err(|source| InitError::CreateLogFile {
            path: path.to_path_buf(),
            source,
        })?;

        let file_layer = tracing_subscriber::fmt::layer()
            .with_writer(Mutex::new(file))
            .with_ansi(false)
            .event_format(CliLogFormatter::new(clock))
            .with_filter(level_filter(verbose));
        // Mirror errors to stderr in `--log` mode, but keep the same
        // terminal-detection ANSI policy as the normal stderr-only logger.
        let stderr_layer = tracing_subscriber::fmt::layer()
            .with_writer(std::io::stderr)
            .with_ansi(stderr_ansi_enabled)
            .event_format(CliLogFormatter::new(clock))
            .with_filter(LevelFilter::ERROR);
        let subscriber = tracing_subscriber::registry()
            .with(file_layer)
            .with(stderr_layer);
        tracing::subscriber::set_global_default(subscriber)?;
    } else {
        let subscriber = tracing_subscriber::fmt()
            .with_max_level(level_filter(verbose))
            .with_writer(std::io::stderr)
            .with_ansi(stderr_ansi_enabled)
            .event_format(CliLogFormatter::new(clock))
            .finish();

        tracing::subscriber::set_global_default(subscriber)?;
    }

    Ok(())
}

fn level_filter(verbose: u8) -> LevelFilter {
    match verbose {
        0 => LevelFilter::WARN,
        1 => LevelFilter::INFO,
        2 => LevelFilter::DEBUG,
        _ => LevelFilter::TRACE,
    }
}

/// A compact formatter tuned for CLI output rather than telemetry output.
#[derive(Debug, Clone)]
struct CliLogFormatter {
    clock: TimestampClock,
}

impl CliLogFormatter {
    fn new(clock: TimestampClock) -> Self {
        Self { clock }
    }

    fn timestamp(&self) -> String {
        self.clock.format_now()
    }

    fn format_prefix(&self, writer: &mut Writer<'_>, level: Level) -> fmt::Result {
        let timestamp = self.timestamp();

        write_styled(writer, level_style(level), level.as_str())?;
        writer.write_char(' ')?;
        write_styled(writer, Style::new().dimmed(), "|")?;
        writer.write_char(' ')?;
        write_styled(writer, Style::new().dimmed(), &timestamp)?;
        writer.write_char(' ')?;
        write_styled(writer, Style::new().dimmed(), "|")?;
        writer.write_char(' ')
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TimestampClock {
    offset: UtcOffset,
    show_utc_suffix: bool,
    #[cfg(test)]
    fixed_timestamp: Option<&'static str>,
}

impl TimestampClock {
    fn detect() -> Self {
        match UtcOffset::current_local_offset() {
            Ok(offset) => Self {
                offset,
                show_utc_suffix: false,
                #[cfg(test)]
                fixed_timestamp: None,
            },
            Err(_) => Self {
                offset: UtcOffset::UTC,
                show_utc_suffix: true,
                #[cfg(test)]
                fixed_timestamp: None,
            },
        }
    }

    #[cfg(test)]
    fn fixed(timestamp: &'static str) -> Self {
        Self {
            offset: UtcOffset::UTC,
            show_utc_suffix: false,
            fixed_timestamp: Some(timestamp),
        }
    }

    fn format_now(self) -> String {
        #[cfg(test)]
        if let Some(timestamp) = self.fixed_timestamp {
            return timestamp.to_string();
        }

        let local_time = OffsetDateTime::now_utc().to_offset(self.offset);
        let mut timestamp = format!(
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            local_time.year(),
            u8::from(local_time.month()),
            local_time.day(),
            local_time.hour(),
            local_time.minute(),
            local_time.second()
        );

        if self.show_utc_suffix {
            timestamp.push_str(" UTC");
        }

        timestamp
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
struct EventBodyVisitor {
    message: Option<String>,
    fields: String,
}

impl EventBodyVisitor {
    fn from_event(event: &Event<'_>) -> Self {
        let mut visitor = Self::default();
        event.record(&mut visitor);
        visitor
    }

    fn record_value(&mut self, field: &Field, value: String) {
        let sanitized_value = if value.contains('\u{1b}') {
            value.replace('\u{1b}', "\\u{1b}")
        } else {
            value
        };

        if field.name() == "message" {
            self.message = Some(sanitized_value);
            return;
        }

        if !self.fields.is_empty() {
            self.fields.push(' ');
        }
        self.fields.push_str(field.name().trim_start_matches("r#"));
        self.fields.push('=');
        self.fields.push_str(&sanitized_value);
    }

    fn write_body(&self, writer: &mut Writer<'_>) -> fmt::Result {
        if let Some(message) = &self.message {
            writer.write_str(message)?;
        }

        if !self.fields.is_empty() {
            if self.message.is_some() {
                writer.write_str("  ")?;
            }
            write_styled(writer, Style::new().dimmed(), &self.fields)?;
        }

        Ok(())
    }
}

impl Visit for EventBodyVisitor {
    fn record_f64(&mut self, field: &Field, value: f64) {
        self.record_value(field, value.to_string());
    }

    fn record_i64(&mut self, field: &Field, value: i64) {
        self.record_value(field, value.to_string());
    }

    fn record_u64(&mut self, field: &Field, value: u64) {
        self.record_value(field, value.to_string());
    }

    fn record_i128(&mut self, field: &Field, value: i128) {
        self.record_value(field, value.to_string());
    }

    fn record_u128(&mut self, field: &Field, value: u128) {
        self.record_value(field, value.to_string());
    }

    fn record_bool(&mut self, field: &Field, value: bool) {
        self.record_value(field, value.to_string());
    }

    fn record_str(&mut self, field: &Field, value: &str) {
        self.record_value(field, value.to_string());
    }

    fn record_error(&mut self, field: &Field, value: &(dyn std::error::Error + 'static)) {
        self.record_value(field, value.to_string());
    }

    fn record_debug(&mut self, field: &Field, value: &dyn fmt::Debug) {
        self.record_value(field, format!("{value:?}"));
    }
}

impl<S, N> FormatEvent<S, N> for CliLogFormatter
where
    S: Subscriber + for<'a> LookupSpan<'a>,
    N: for<'a> FormatFields<'a> + 'static,
{
    fn format_event(
        &self,
        _ctx: &FmtContext<'_, S, N>,
        mut writer: Writer<'_>,
        event: &Event<'_>,
    ) -> fmt::Result {
        let body = EventBodyVisitor::from_event(event);
        self.format_prefix(&mut writer, *event.metadata().level())?;
        body.write_body(&mut writer)?;
        writeln!(writer)
    }
}

fn level_style(level: Level) -> Style {
    match level {
        Level::ERROR => Color::Red.bold(),
        Level::WARN => Color::Yellow.bold(),
        Level::INFO => Color::Cyan.bold(),
        Level::DEBUG => Color::Blue.bold(),
        Level::TRACE => Color::Purple.bold(),
    }
}

fn write_styled(writer: &mut Writer<'_>, style: Style, text: &str) -> fmt::Result {
    if writer.has_ansi_escapes() {
        write!(writer, "{}", style.paint(text))
    } else {
        writer.write_str(text)
    }
}

#[cfg(test)]
mod tests {
    use super::{CliLogFormatter, TimestampClock, level_filter};
    use std::io;
    use std::sync::{Arc, Mutex};
    use tracing_subscriber::filter::LevelFilter;
    use tracing_subscriber::layer::SubscriberExt;

    const FIXED_TIMESTAMP: &str = "2026-04-23 11:03:42";

    struct BufferWriter(Arc<Mutex<Vec<u8>>>);

    impl io::Write for BufferWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().expect("buffer mutex poisoned").extend(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn capture_log(ansi_enabled: bool, emit: impl FnOnce()) -> String {
        let buffer = Arc::new(Mutex::new(Vec::new()));
        let writer_buffer = Arc::clone(&buffer);
        let layer = tracing_subscriber::fmt::layer()
            .with_writer(move || BufferWriter(Arc::clone(&writer_buffer)))
            .with_ansi(ansi_enabled)
            .event_format(CliLogFormatter::new(TimestampClock::fixed(FIXED_TIMESTAMP)));
        let subscriber = tracing_subscriber::registry().with(layer);

        tracing::subscriber::with_default(subscriber, emit);

        let bytes = buffer.lock().expect("buffer mutex poisoned").clone();
        String::from_utf8(bytes).expect("log output should be valid UTF-8")
    }

    #[test]
    fn verbose_levels_enable_warn_by_default_and_trace_at_three_v() {
        assert_eq!(level_filter(0), LevelFilter::WARN);
        assert_eq!(level_filter(1), LevelFilter::INFO);
        assert_eq!(level_filter(2), LevelFilter::DEBUG);
        assert_eq!(level_filter(3), LevelFilter::TRACE);
        assert_eq!(level_filter(4), LevelFilter::TRACE);
        assert_eq!(level_filter(u8::MAX), LevelFilter::TRACE);
    }

    #[test]
    fn ansi_line_contains_escape_sequences() {
        let line = capture_log(true, || tracing::error!("Command failed"));

        assert!(line.contains('\u{1b}'));
    }

    #[test]
    fn plain_line_has_no_escape_sequences() {
        let line = capture_log(false, || tracing::error!("Command failed"));

        assert!(!line.contains('\u{1b}'));
    }

    #[test]
    fn message_only_event_uses_level_timestamp_message_layout() {
        let line = capture_log(false, || tracing::info!("Loaded input alignment"));

        assert_eq!(
            line,
            "INFO | 2026-04-23 11:03:42 | Loaded input alignment\n"
        );
    }

    #[test]
    fn message_and_fields_keep_message_before_flat_fields() {
        let line = capture_log(false, || {
            tracing::trace!(
                iteration = 3_u64,
                working_set_count = 18_u64,
                "Selected heuristic candidate"
            );
        });

        assert_eq!(
            line,
            "TRACE | 2026-04-23 11:03:42 | Selected heuristic candidate  iteration=3 working_set_count=18\n"
        );
    }

    #[test]
    fn field_only_event_has_no_extra_body_padding() {
        let line = capture_log(false, || tracing::warn!(missing = "seq1"));

        assert_eq!(line, "WARN | 2026-04-23 11:03:42 | missing=seq1\n");
    }

    #[test]
    fn raw_identifier_prefixes_are_removed_from_field_names() {
        let line = capture_log(false, || tracing::info!(r#type = "dna", "Raw field"));

        assert_eq!(line, "INFO | 2026-04-23 11:03:42 | Raw field  type=dna\n");
    }

    #[test]
    fn ansi_escape_characters_are_escaped_in_message_and_fields() {
        let line = capture_log(false, || {
            tracing::warn!(name = "bad\u{1b}[31m", "Hello\u{1b}[0m");
        });

        assert_eq!(
            line,
            "WARN | 2026-04-23 11:03:42 | Hello\\u{1b}[0m  name=bad\\u{1b}[31m\n"
        );
    }
}
