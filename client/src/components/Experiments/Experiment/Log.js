import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "@material-ui/core/styles";

class Log extends Component {
  render() {
    return (
      <Container>
        {this.renderEntry(
          "create",
          this.props.created,
          this.statuses("info"),
          "Experiment created"
        )}
        {Object.keys(this.props.pipeline).map(this.renderAction.bind(this))}
        {this.props.done &&
          this.renderEntry(
            "done",
            this.props.done,
            this.statuses("info"),
            "Experiment completed"
          )}
      </Container>
    );
  }

  statuses(key) {
    const { primary, grey, error, warning } = this.props.theme.palette;
    const statuses = {
      completed: { text: "DONE", color: primary.main },
      running: { text: "RUNNING", color: warning.main },
      waiting: { text: "WAITING", color: grey[500] },
      info: { text: "INFO", color: grey[500] },
      error: { text: "ERROR", color: error.main },
      interrupted: { text: "INTERRUPTED", color: warning.main }
    };
    return statuses[key];
  }

  renderEntry(key, time, status, content) {
    return (
      <Entry key={`entry-${key}`}>
        {`[${formatTime(time)}] `}
        <Status color={status.color}>{status.text}</Status>
        {` ${content}`}
      </Entry>
    );
  }

  renderAction(actionName) {
    const action = this.props.pipeline[actionName];
    return this.renderEntry(
      actionName,
      this.actionTime(action),
      this.actionStatus(action),
      this.formatActionContent(action, actionName)
    );
  }

  actionTime(action) {
    return action.completed
      ? action.completed
      : action.error
        ? action.error
        : action.interrupted
          ? action.interrupted
          : action.started
            ? action.started
            : null;
  }

  actionStatus(action) {
    return action.error
      ? this.statuses("error")
      : action.completed
        ? this.statuses("completed")
        : action.started
          ? action.interrupted
            ? this.statuses("interrupted")
            : this.statuses("running")
          : this.statuses("waiting");
  }

  formatActionContent(action, actionName) {
    const upperCaseName =
      actionName.substr(0, 1).toUpperCase() + actionName.substr(1);
    return action.error
      ? this.props.error
      : action.cached
        ? `Load ${actionName} from disk`
        : upperCaseName;
  }
}

const formatTime = time => {
  const noTime = ["----", "--", "--", "--", "--", "--"];
  time = time || noTime;
  return `${time[0]}/${pad(time[1])}/${pad(time[2])} ${pad(time[3])}:${pad(
    time[4]
  )}:${pad(time[5])}`;
};

const pad = number => {
  number = number.toString();
  return number.length === 2 ? number : "0" + number;
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
