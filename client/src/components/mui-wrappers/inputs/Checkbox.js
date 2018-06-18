import React, { Component } from "react";
import styled from "styled-components";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import Checkbox from "@material-ui/core/Checkbox";

export default class extends Component {
  render() {
    return (
      <StyledFormControlLabel
        control={<Checkbox {...this.props} color="primary" />}
        label={this.props.label}
      />
    );
  }
}

const StyledFormControlLabel = styled(FormControlLabel)`
  white-space: nowrap;
`;
