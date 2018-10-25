import React, { Component } from "react";
import styled from "styled-components";
import TextField from "@material-ui/core/TextField";
import { width, marginTop, marginRight, marginBottom } from "./constants";

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
        width={this.props.width || width}
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
  margin-right: ${marginRight}px !important;
  margin-top: ${marginTop}px !important;
  margin-bottom: ${marginBottom}px !important;
`;
