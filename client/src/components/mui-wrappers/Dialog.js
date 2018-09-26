import React, { Component } from "react";
import styled from "styled-components";
import Dialog from "@material-ui/core/Dialog";
import DialogActions from "@material-ui/core/DialogActions";
import DialogTitle from "@material-ui/core/DialogTitle";
import DialogContent from "@material-ui/core/DialogContent";
import Button from "@material-ui/core/Button";

export default class extends Component {
  render() {
    return (
      <Dialog open={this.props.open} aria-labelledby="dialog-title">
        <DialogTitle id="dialog-title">{this.props.title}</DialogTitle>
        <StyledDialogContent>{this.props.children}</StyledDialogContent>
        {this.props.actions && (
          <DialogActions>
            {this.props.actions.map(this.renderAction)}
          </DialogActions>
        )}
      </Dialog>
    );
  }

  renderAction(action, index) {
    return (
      <Button
        key={`action-${index}`}
        onClick={action.onClick}
        color={action.color || "default"}
        disabled={action.disabled}
      >
        {action.name}
      </Button>
    );
  }
}

const StyledDialogContent = styled(DialogContent)`
  display: flex;
  align-items: center;
`;
