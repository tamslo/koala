import React, { Component } from "react";
import styled from "styled-components";
import TextField from "@material-ui/core/TextField";

export default class extends Component {
  render() {
    return (
      <FixedWidthTextField
        label={this.props.label}
        value={this.props.value}
        onChange={this.props.onChange}
        margin="normal"
        select={this.props.select}
        fullWidth={this.props.fullWidth}
        width={this.props.width || 200}
        type={this.props.type || "text"}
        disabled={this.props.disabled}
      >
        {this.props.children}
      </FixedWidthTextField>
    );
  }
}

const FixedWidthTextField = styled(TextField)`
  width: ${props => props.width}px;
  margin-right: 20px !important;
  .MuiInput-input-134,
  .jss-134 {
    display: flex !important;
    justify-content: space-between !important;
  }
`;
