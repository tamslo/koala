import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "material-ui/styles";

class Experiment extends Component {
  render() {
    return (
      <div>
        <Entry>{`Data URL: ${this.props.dataset}`}</Entry>
        <Entry>{`Aligner: ${this.props.aligner}`}</Entry>
        {this.renderLog()}
      </div>
    );
  }

  renderLog() {
    return <Log>{this.props.log.map(this.renderLogEntry.bind(this))}</Log>;
  }

  renderLogEntry(entry, index) {
    const { theme } = this.props;
    const { primary, error, warning, text } = theme.palette;
    const texts = {
      create: "Create experiment",
      cached: "Load data from cache",
      dataset: "Download data",
      update: "Update experiment"
    };

    const time = entry.completed
      ? formatTime(entry.completed)
      : entry.error ? formatTime(entry.error) : formatTime(entry.started);

    const status =
      entry.action === "interrupted" || entry.action === "done"
        ? { text: "INFO", color: warning.main }
        : entry.completed
          ? { text: "DONE", color: primary.main }
          : entry.error
            ? { text: "ERROR", color: error.main }
            : { text: "START", color: text.secondary };

    const content = entry.error
      ? this.props.error
      : texts[entry.action] || capitalize(entry.action);

    return (
      <Entry key={`log-entry-${index}`}>
        {`[${time}] [`}
        <Status color={status.color}>{status.text}</Status>
        {`] ${content}`}
      </Entry>
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

const Log = styled.div`
  padding-top: 12px;
  font-family: monospace;
`;

const Entry = styled.div`
  margin-bottom: 12px;
`;

const Status = styled.span`
  color: ${props => props.color};
`;

export default withTheme()(Experiment);
