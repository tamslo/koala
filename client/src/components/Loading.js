import React, { Component } from "react";
import styled from "styled-components";
import Dialog from "@material-ui/core/Dialog";
import DialogActions from "@material-ui/core/DialogActions";
import DialogTitle from "@material-ui/core/DialogTitle";
import DialogContent from "@material-ui/core/DialogContent";
import DialogContentText from "@material-ui/core/DialogContentText";
import CircularProgress from "@material-ui/core/CircularProgress";
import Button from "@material-ui/core/Button";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = { open: true };
  }

  handleClose() {
    this.setState({ open: false });
  }

  render() {
    const title = this.props.error ? "Something went wrong" : "Stay with us";
    return (
      <Dialog
        open={this.state.open}
        aria-labelledby="dialog-title"
        aria-describedby="alert-dialog-description"
      >
        <DialogTitle id="dialog-title">{title}</DialogTitle>
        <StyledDialogContent>
          <DialogContentText id="alert-dialog-description">
            {this.props.content}
          </DialogContentText>
          {!this.props.error && <StyledCircularProgress />}
        </StyledDialogContent>
        {this.props.error && (
          <DialogActions>
            <Button onClick={this.handleClose.bind(this)} color="default">
              Cancel
            </Button>
            <Button onClick={this.props.retry} color="primary">
              Retry
            </Button>
          </DialogActions>
        )}
      </Dialog>
    );
  }
}

const StyledDialogContent = styled(DialogContent)`
  display: flex;
  align-items: center;
`;

const StyledCircularProgress = styled(CircularProgress)`
  margin-left: 10px;
  height: 30px !important;
  width: 30px !important;
`;
