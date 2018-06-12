import React, { Component } from "react";
import TextField from "./Text";

export default class extends Component {
  render() {
    return (
      <TextField {...this.props} type="number">
        {this.props.children}
      </TextField>
    );
  }
}
