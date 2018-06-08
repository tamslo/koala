import React, { Component } from "react";
import styled from "styled-components";
import DialogContentText from "@material-ui/core/DialogContentText";
import CircularProgress from "@material-ui/core/CircularProgress";
import Dialog from "./mui-wrappers/Dialog";

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
    const actions = [
      {
        name: "Cancel",
        onClick: this.props.retry,
        disabled: !this.props.error,
        color: "primary"
      },
      {
        name: "Retry",
        onClick: this.handleClose.bind(this),
        disabled: !this.props.error
      }
    ];
    return (
      <Dialog open={this.state.open} title={title} actions={actions}>
        <DialogContentText id="alert-dialog-description">
          {this.props.content}
        </DialogContentText>
        {!this.props.error && <StyledCircularProgress />}
      </Dialog>
    );
  }
}

const StyledCircularProgress = styled(CircularProgress)`
  margin-left: 10px;
  height: 30px !important;
  width: 30px !important;
`;
