import React, { Component } from "react";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import Checkbox from "@material-ui/core/Checkbox";

export default class extends Component {
  render() {
    return (
      <FormControlLabel
        control={
          <Checkbox
            checked={this.props.checked}
            onChange={this.props.onChange}
            color="primary"
          />
        }
        label={this.props.label}
      />
    );
  }
}
