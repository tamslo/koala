import React, { Component } from "react";
import List, {
  ListItem,
  ListItemSecondaryAction,
  ListItemText
} from "material-ui/List";
import IconButton from "material-ui/IconButton";
import DeleteIcon from "@material-ui/icons/Delete";
import Card from "./Card";

export default class extends Component {
  render() {
    console.log(this.props.experiments);
    const { experiments } = this.props;
    return (
      <Card title="Executed experiments">
        <List>
          {Object.keys(experiments).map((experimentId, index) => (
            <ListItem key={experimentId} button>
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
      </Card>
    );
  }

  primaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return experiment["name"];
  }

  secondaryText(experimentId) {
    const experiment = this.props.experiments[experimentId];
    return experiment["created"];
  }
}
