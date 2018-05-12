import React, { Component } from "react";
import styled from "styled-components";
import TextField from "material-ui/TextField";
import MenuItem from "material-ui/Menu/MenuItem";
import Button from "material-ui/Button";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = {
      name: "New Experiment",
      dataset: "https://www.example.com",
      aligner: "star"
    };
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  render() {
    const { aligners } = this.props;
    return (
      <Container>
        <FixedWidthTextField
          label="Name"
          value={this.state.name}
          onChange={this.handleChange("name")}
          margin="normal"
        />
        <DataInput
          label="Data URL"
          value={this.state.dataset}
          onChange={this.handleChange("dataset")}
          margin="normal"
        />
        <FixedWidthTextField
          select
          label="Aligner"
          value={this.state.aligner || ""}
          onChange={this.handleChange("aligner")}
          margin="normal"
        >
          <MenuItem value="">
            <em>None</em>
          </MenuItem>
          {Object.keys(aligners).map(aligner => (
            <MenuItem key={aligner} value={aligner}>
              {aligners[aligner].name}
            </MenuItem>
          ))}
        </FixedWidthTextField>
        <Button
          color="primary"
          onClick={() => this.props.addExperiment(this.state)}
          disabled={!this.canRun()}
          size="large"
        >
          Add
        </Button>
      </Container>
    );
  }

  canRun() {
    const nameValid = this.state.name !== "";
    const dataValid = this.state.dataset !== "";
    const alignerValid = Object.keys(this.props.aligners).includes(
      this.state.aligner
    );
    return nameValid && dataValid && alignerValid;
  }
}

const Container = styled.div`
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  padding-left: 12px;
`;

const DataInput = styled(TextField)`
  flex-grow: 1;
  min-width: 200px;
  margin-right: 20px !important;
`;

const FixedWidthTextField = styled(TextField)`
  width: 200px;
  margin-right: 20px !important;
`;
