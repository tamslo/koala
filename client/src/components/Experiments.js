import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "material-ui/styles";
import List, {
  ListItem,
  ListItemSecondaryAction,
  ListItemText
} from "material-ui/List";
import IconButton from "material-ui/IconButton";
import Collapse from "material-ui/transitions/Collapse";
import DeleteIcon from "@material-ui/icons/Delete";
import ExpandLess from "@material-ui/icons/ExpandLess";
import ExpandMore from "@material-ui/icons/ExpandMore";
import Card from "./Card";
import Experiment from "./Experiment";

class Experiments extends Component {
  constructor(props) {
    super(props);
    this.state = { open: [] };
  }

  toggleExperiment(id) {
    if (this.state.open.includes(id)) {
      this.setState({
        open: this.state.open.filter(experimentId => experimentId !== id)
      });
    } else {
      this.setState({ open: [...this.state.open, id] });
    }
  }

  render() {
    return (
      <Card title="Executed experiments">
        <List>
          {Object.keys(this.props.experiments)
            .reverse()
            .map(this.renderExperiment.bind(this))}
        </List>
      </Card>
    );
  }

  renderExperiment(experimentId) {
    const { experiments } = this.props;
    const experiment = experiments[experimentId];
    return (
      <div key={experiment.id}>
        {this.renderListItem(experiment)}
        <Collapse
          in={this.state.open.includes(experiment.id)}
          timeout="auto"
          unmountOnExit
        >
          <InsetListItem>
            <Experiment {...experiment} />
          </InsetListItem>
        </Collapse>
      </div>
    );
  }

  renderListItem(experiment) {
    const { deleteExperiment } = this.props;
    return (
      <ListItem
        key={experiment.id}
        button
        onClick={() => {
          this.toggleExperiment(experiment.id);
        }}
      >
        <ListItemText
          primary={this.primaryText(experiment.id)}
          secondary={this.secondaryText(experiment.id)}
        />
        <ListItemSecondaryAction>
          <IconButton
            aria-label="Delete"
            onClick={() => {
              deleteExperiment(experiment.id);
            }}
            disabled={
              !(experiment.done || experiment.interrupted || experiment.error)
            }
          >
            <DeleteIcon />
          </IconButton>
          <IconButton
            aria-label="Collapse"
            onClick={() => {
              this.toggleExperiment(experiment.id);
            }}
          >
            {this.state.open.includes(experiment.id) ? (
              <ExpandLess />
            ) : (
              <ExpandMore />
            )}
          </IconButton>
        </ListItemSecondaryAction>
      </ListItem>
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
          : "Error"}
      </Status>
    );
  }
}

const Status = styled.span`
  color: ${props => props.color};
`;

const InsetListItem = styled(ListItem)`
  padding-left: 32px !important;
`;

export default withTheme()(Experiments);
