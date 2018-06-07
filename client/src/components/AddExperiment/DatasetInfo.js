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
      <Dialog
        open={this.props.open}
        aria-labelledby="dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="dialog-title">{this.props.dataset.name}</DialogTitle>
        <StyledDialogContent>
          <div>{`URL: ${this.props.dataset.url}`}</div>
        </StyledDialogContent>
        <DialogActions>
          <Button onClick={() => this.props.close} color="default" disabled>
            Delete
          </Button>
          <Button onClick={this.props.close} color="default">
            Close
          </Button>
        </DialogActions>
      </Dialog>
    );
  }

  addDataset() {
    this.setState(this.initialState(), this.props.addDataset(this.state));
  }
}

const StyledDialogContent = styled(DialogContent)`
  display: flex;
  align-items: center;
`;
