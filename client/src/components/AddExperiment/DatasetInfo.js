import React, { Component } from "react";
// import styled from "styled-components";
import Dialog from "../mui-wrappers/Dialog";

export default class extends Component {
  render() {
    const actions = [
      {
        name: "Delete",
        onClick: () => this.props.close,
        disabled: true
      },
      {
        name: "Close",
        onClick: this.props.close,
        color: "primary"
      }
    ];
    return (
      <Dialog
        open={this.props.open}
        title={this.props.dataset.name}
        actions={actions}
      >
        <div>{`URL: ${this.props.dataset.url}`}</div>
      </Dialog>
    );
  }
}
