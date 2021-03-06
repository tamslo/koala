import React, { Component } from "react";
import styled from "styled-components";
import { withTheme } from "@material-ui/core/styles";
import ListSubheader from "@material-ui/core/ListSubheader";
import List from "@material-ui/core/List";
import ListItem from "@material-ui/core/ListItem";
import ListItemSecondaryAction from "@material-ui/core/ListItemSecondaryAction";
import ListItemText from "@material-ui/core/ListItemText";
import IconButton from "@material-ui/core/IconButton";
import Collapse from "@material-ui/core/Collapse";
import DeleteIcon from "@material-ui/icons/Delete";
import RefreshIcon from "@material-ui/icons/Refresh";
import ExpandLessIcon from "@material-ui/icons/ExpandLess";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import Experiment from "./Experiment";
import { SERVER_URL } from "../../api";
import constants from "../../constants.json";

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
      <List
        subheader={<StyledHeader component="div">Experiments</StyledHeader>}
      >
        {Object.keys(this.props.experiments)
          .reverse()
          .map(this.renderExperiment.bind(this))}
      </List>
    );
  }

  renderExperiment(experimentId) {
    const { experiments, services, datasets, references } = this.props;
    const experiment = experiments[experimentId];
    const experimentProps = {
      services,
      datasets,
      references,
      SERVER_URL,
      ...experiment
    };
    return (
      <div key={experiment.id}>
        {this.renderListItem(experiment)}
        <Collapse
          in={this.state.open.includes(experiment.id)}
          timeout="auto"
          unmountOnExit
        >
          <InsetListItem>
            <Experiment {...experimentProps} />
          </InsetListItem>
        </Collapse>
      </div>
    );
  }

  renderListItem(experiment) {
    const { deleteExperiment, retryExperiment } = this.props;
    const retryEnabled =
      (constants.retryDoneExperiment ||
        experiment.status === constants.experiment.ERROR) &&
      experiment.status !== constants.experiment.WAITING &&
      experiment.status !== constants.experiment.RUNNING;
    return (
      <StyledListItem
        key={experiment.id}
        button
        onClick={() => {
          this.toggleExperiment(experiment.id);
        }}
        disableGutters
      >
        <ListItemText
          primary={this.primaryText(experiment)}
          secondary={this.secondaryText(experiment)}
        />
        <ListItemSecondaryAction>
          {retryEnabled && (
            <IconButton
              aria-label="Retry"
              onClick={() => retryExperiment(experiment)}
            >
              <RefreshIcon />
            </IconButton>
          )}
          {experiment.status === constants.experiment.DONE && (
            <IconButton
              aria-label="Download"
              href={SERVER_URL + "/export?experiment=" + experiment.id}
            >
              <DownloadIcon />
            </IconButton>
          )}
          <IconButton
            aria-label="Delete"
            onClick={() => deleteExperiment(experiment.id)}
            disabled={experiment.status === constants.experiment.RUNNING}
          >
            <DeleteIcon />
          </IconButton>
          <IconButton
            aria-label="Collapse"
            onClick={() => this.toggleExperiment(experiment.id)}
          >
            {this.state.open.includes(experiment.id) ? (
              <ExpandLessIcon />
            ) : (
              <ExpandMoreIcon />
            )}
          </IconButton>
        </ListItemSecondaryAction>
      </StyledListItem>
    );
  }

  primaryText(experiment) {
    return experiment["name"];
  }

  statusColor(experiment) {
    const { primary, error } = this.props.theme.palette;
    const { status } = experiment;
    return status === constants.experiment.ERROR
      ? error.main
      : status === constants.experiment.RUNNING
        ? primary.main
        : "inherit";
  }

  secondaryText(experiment) {
    return (
      <Status color={this.statusColor(experiment)}>
        {capitalize(experiment.status)}
      </Status>
    );
  }
}

const capitalize = string => string.charAt(0).toUpperCase() + string.substr(1);

const Status = styled.span`
  color: ${props => props.color};
`;

const StyledHeader = styled(ListSubheader)`
  padding-left: 12px !important;
`;

const StyledListItem = styled(ListItem)`
  padding-left: 12px !important;
`;

const InsetListItem = styled(ListItem)`
  padding-left: 32px !important;
`;

export default withTheme()(Experiments);
