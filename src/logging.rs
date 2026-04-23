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
use tracing_subscriber::filter::filter_fn;
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
            .event_format(CliLogFormatter::new(clock, false))
            .with_filter(level_filter(verbose));
        // Mirror errors to stderr in `--log` mode, but keep the same
        // terminal-detection ANSI policy as the normal stderr-only logger.
        let stderr_layer = tracing_subscriber::fmt::layer()
            .with_writer(std::io::stderr)
            .with_ansi(stderr_ansi_enabled)
            .event_format(CliLogFormatter::new(clock, stderr_ansi_enabled))
            .with_filter(filter_fn(|metadata| metadata.level() == &Level::ERROR));
        let subscriber = tracing_subscriber::registry()
            .with(file_layer)
            .with(stderr_layer);
        tracing::subscriber::set_global_default(subscriber)?;
    } else {
        let subscriber = tracing_subscriber::fmt()
            .with_max_level(level_filter(verbose))
            .with_writer(std::io::stderr)
            .with_ansi(stderr_ansi_enabled)
            .event_format(CliLogFormatter::new(clock, stderr_ansi_enabled))
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
    ansi_enabled: bool,
}

impl CliLogFormatter {
    fn new(clock: TimestampClock, ansi_enabled: bool) -> Self {
        Self {
            clock,
            ansi_enabled,
        }
    }

    fn timestamp(&self) -> String {
        self.clock.format_now()
    }

    fn format_line(&self, level: Level, message: Option<&str>, fields: Option<&str>) -> String {
        let timestamp = self.timestamp();
        let mut rendered = String::new();

        rendered.push_str(&styled_text(
            level_style(level),
            &level.to_string(),
            self.ansi_enabled,
        ));
        rendered.push(' ');
        rendered.push_str(&styled_text(Style::new().dimmed(), "|", self.ansi_enabled));
        rendered.push(' ');
        rendered.push_str(&styled_text(
            Style::new().dimmed(),
            &timestamp,
            self.ansi_enabled,
        ));
        rendered.push(' ');
        rendered.push_str(&styled_text(Style::new().dimmed(), "|", self.ansi_enabled));
        rendered.push(' ');

        if let Some(message) = message {
            rendered.push_str(message);
        }

        if let Some(fields) = fields.filter(|fields| !fields.is_empty()) {
            if message.is_some() {
                rendered.push_str("  ");
            }
            rendered.push_str(&styled_text(
                Style::new().dimmed(),
                fields,
                self.ansi_enabled,
            ));
        }

        rendered
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TimestampClock {
    offset: UtcOffset,
    show_utc_suffix: bool,
}

impl TimestampClock {
    fn detect() -> Self {
        match UtcOffset::current_local_offset() {
            Ok(offset) => Self {
                offset,
                show_utc_suffix: false,
            },
            Err(_) => Self {
                offset: UtcOffset::UTC,
                show_utc_suffix: true,
            },
        }
    }

    fn format_now(self) -> String {
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
struct EventContent {
    message: Option<String>,
    fields: Vec<EventField>,
}

impl EventContent {
    fn from_event(event: &Event<'_>) -> Self {
        let mut visitor = EventContentVisitor::default();
        event.record(&mut visitor);
        visitor.content
    }

    fn fields_text(&self) -> String {
        let mut rendered = String::new();

        for (index, field) in self.fields.iter().enumerate() {
            if index > 0 {
                rendered.push(' ');
            }
            rendered.push_str(&field.name);
            rendered.push('=');
            rendered.push_str(&field.value);
        }

        rendered
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct EventField {
    name: String,
    value: String,
}

#[derive(Debug, Default)]
struct EventContentVisitor {
    content: EventContent,
}

impl EventContentVisitor {
    fn record_value(&mut self, field: &Field, value: String) {
        let sanitized_value = value.replace('\u{1b}', "\\u{1b}");

        if field.name() == "message" {
            self.content.message = Some(sanitized_value);
            return;
        }

        self.content.fields.push(EventField {
            name: field.name().trim_start_matches("r#").to_string(),
            value: sanitized_value,
        });
    }
}

impl Visit for EventContentVisitor {
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
        let content = EventContent::from_event(event);
        let fields = (!content.fields.is_empty()).then(|| content.fields_text());
        let line = self.format_line(
            *event.metadata().level(),
            content.message.as_deref(),
            fields.as_deref(),
        );

        writer.write_str(&line)?;
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

fn styled_text(style: Style, text: &str, ansi_enabled: bool) -> String {
    if ansi_enabled {
        style.paint(text).to_string()
    } else {
        text.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::{CliLogFormatter, TimestampClock, level_filter};
    use time::UtcOffset;
    use tracing::Level;
    use tracing_subscriber::filter::LevelFilter;

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
        let formatter = CliLogFormatter::new(
            TimestampClock {
                offset: UtcOffset::UTC,
                show_utc_suffix: false,
            },
            true,
        );
        let line = formatter.format_line(Level::ERROR, Some("Command failed"), None);

        assert!(line.contains('\u{1b}'));
    }

    #[test]
    fn plain_line_has_no_escape_sequences() {
        let formatter = CliLogFormatter::new(
            TimestampClock {
                offset: UtcOffset::UTC,
                show_utc_suffix: false,
            },
            false,
        );
        let line = formatter.format_line(Level::ERROR, Some("Command failed"), None);

        assert!(!line.contains('\u{1b}'));
    }
}
