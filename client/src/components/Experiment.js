import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "@material-ui/core/styles";
import IconButton from "@material-ui/core/IconButton";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import { SERVER_URL } from "../api";

class Experiment extends Component {
  render() {
    return (
      <div>
        <Entry>
          {`Data Set: ${this.props.datasets[this.props.dataset].name}`}
          {this.renderDownloadButton("dataset")}
        </Entry>
        <Entry>
          {`Aligner: ${
            this.props.services.find(
              service => service.id === this.props.alignment
            ).name
          }`}
          {this.renderDownloadButton("alignment")}
        </Entry>
        {this.renderLog()}
      </div>
    );
  }

  renderDownloadButton(key) {
    const path = this.props.files[key];
    return path ? (
      <StyledIconButton
        aria-label={"Download"}
        href={SERVER_URL + "/download?path=" + path}
      >
        <DownloadIcon />
      </StyledIconButton>
    ) : null;
  }

  renderLog() {
    return <Log>{this.props.log.map(this.renderLogEntry.bind(this))}</Log>;
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
    return capitalize(string);
  }
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

const StyledIconButton = styled(IconButton)`
  margin-left: 12px !important;
`;

export default withTheme()(Experiment);
