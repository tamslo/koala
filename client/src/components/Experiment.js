import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "material-ui/styles";

class Experiment extends Component {
  render() {
    return <Container>{this.renderLog()}</Container>;
  }

  renderLog() {
    return (
      <div>
        {Object.keys(this.props.log).map(this.renderLogEntry.bind(this))}
      </div>
    );
  }

  renderLogEntry(action, index) {
    const { log, theme } = this.props;
    const { primary, error, warning, text } = theme.palette;
    const entry = log[action];
    const texts = {
      create: "Create experiment",
      cached: "Load data from cache",
      dataset: "Download data"
    };

    const time = entry.completed
      ? formatTime(entry.completed)
      : entry.error ? formatTime(entry.error) : formatTime(entry.started);

    const status =
      action === "interrupted"
        ? { text: "INFO", color: warning.main }
        : entry.completed
          ? { text: "DONE", color: primary.main }
          : entry.error
            ? { text: "ERROR", color: error.main }
            : { text: "START", color: text.secondary };

    return (
      <LogEntry key={`log-entry-${index}`} index={index}>
        {`[${time}] [`}
        <Status color={status.color}>{status.text}</Status>
        {`] ${this.props.error || texts[action] || capitalize(action)}`}
      </LogEntry>
    );
  }
}

const formatTime = time => {
  return `${time[0]}/${pad(time[1])}/${pad(time[2])} ${pad(time[3])}:${pad(
    time[4]
  )}:${pad(time[5])}`;
};

const pad = number => {
  number = number.toString();
  return number.length === 2 ? number : "0" + number;
};

const capitalize = string => {
  return string.substr(0, 1).toUpperCase() + string.substr(1);
};

const Container = styled.div`
  padding: 12px;
`;

const LogEntry = styled.div`
  margin-top: ${props => (props.index === 0 ? 0 : 12)}px;
  font-family: monospace;
`;

const Status = styled.span`
  color: ${props => props.color};
`;

export default withTheme()(Experiment);
