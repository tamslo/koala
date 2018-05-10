import React, { Component } from "react";
import List, {
  ListItem,
  ListItemSecondaryAction,
  ListItemText
} from "material-ui/List";
import Divider from "material-ui/Divider";
import IconButton from "material-ui/IconButton";
import DeleteIcon from "@material-ui/icons/Delete";

export default class extends Component {
  render() {
    console.log(this.props.experiments);
    const { experiments } = this.props;
    return (
      <List>
        {Object.keys(experiments).map((experimentId, index) => (
          <ListItem key={experimentId}>
            <ListItemText
              primary={this.primaryText(experimentId)}
              secondary={this.secondaryText(experimentId)}
            />
            <ListItemSecondaryAction>
              <IconButton aria-label="Delete">
                <DeleteIcon />
              </IconButton>
            </ListItemSecondaryAction>
          </ListItem>
        ))}
      </List>
    );
  }

  primaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return experiment["dataset"];
  }

  secondaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return experiment["created"];
  }
}
