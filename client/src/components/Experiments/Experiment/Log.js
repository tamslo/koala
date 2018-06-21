import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "@material-ui/core/styles";

class Log extends Component {
  render() {
    return (
      <Container>
        {this.props.log.map(this.renderLogEntry.bind(this))}
      </Container>
    );
  }

  renderLogEntry(entry, index) {
    const { theme } = this.props;
    const { primary, error, warning, text } = theme.palette;
    const texts = {
      create: "Create experiment",
      dataset: "Download data",
      update: "Update experiment"
    };

    const time = entry.completed
      ? formatTime(entry.completed)
      : entry.error
        ? formatTime(entry.error)
        : formatTime(entry.started);

    const status =
      entry.action === "done"
        ? { text: "INFO", color: text.secondary }
        : entry.error
          ? { text: "ERROR", color: error.main }
          : entry.interrupted
            ? { text: "INTERRUPT", color: warning.main }
            : { text: "DONE", color: primary.main };

    const content = entry.error
      ? this.props.error
      : texts[entry.action] || formatAction(entry.action);

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

const formatAction = string => {
  if (string.endsWith("_loaded")) {
    return `Load ${string.split("_")[0]} from disk`;
  } else {
    return string.substr(0, 1).toUpperCase() + string.substr(1);
  }
};

const Container = styled.div`
  padding-top: 12px;
  font-family: monospace;
`;

const Entry = styled.div`
  margin-bottom: 10px;
`;

const Status = styled.span`
  color: ${props => props.color};
`;

export default withTheme()(Log);
