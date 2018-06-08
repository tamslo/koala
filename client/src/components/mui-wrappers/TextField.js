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
      >
        {this.props.children}
      </FixedWidthTextField>
    );
  }
}

const FixedWidthTextField = styled(TextField)`
  width: 200px;
  margin-right: 20px !important;
`;
