import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "material-ui/styles";
import List, {
  ListItem,
  ListItemSecondaryAction,
  ListItemText
} from "material-ui/List";
import IconButton from "material-ui/IconButton";
import DeleteIcon from "@material-ui/icons/Delete";
import Card from "./Card";

class Experiments extends Component {
  render() {
    const { experiments, deleteExperiment } = this.props;
    return (
      <Card title="Executed experiments">
        <List>
          {Object.keys(experiments)
            .reverse()
            .map((experimentId, index) => (
              <ListItem key={experimentId} button>
                <ListItemText
                  primary={this.primaryText(experimentId)}
                  secondary={this.secondaryText(experimentId)}
                />
                <ListItemSecondaryAction>
                  <IconButton
                    aria-label="Delete"
                    onClick={() => {
                      deleteExperiment(experimentId);
                    }}
                    disabled={
                      !(
                        experiments[experimentId].done ||
                        experiments[experimentId].interrupted ||
                        experiments[experimentId].error
                      )
                    }
                  >
                    <DeleteIcon />
                  </IconButton>
                </ListItemSecondaryAction>
              </ListItem>
            ))}
        </List>
      </Card>
    );
  }

  primaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return experiment["name"];
  }

  statusColor(experimentId) {
    const experiment = this.props.experiments[experimentId];
    const { primary, error, warning } = this.props.theme.palette;
    return experiment.done
      ? "inherit"
      : experiment.interrupted
        ? warning.main
        : experiment.error ? error.main : primary.main;
  }

  secondaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return (
      <Status color={this.statusColor(experimentId)}>
        {!experiment.error
          ? experiment.done
            ? "Done"
            : experiment.interrupted ? "Interrupted" : "Running"
          : experiment.error}
      </Status>
    );
  }
}

const Status = styled.span`
  color: ${props => props.color};
`;

export default withTheme()(Experiments);
