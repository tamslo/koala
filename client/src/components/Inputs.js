import React, { Component } from "react";
import styled from "styled-components";
import TextField from "material-ui/TextField";
import MenuItem from "material-ui/Menu/MenuItem";
import Button from "material-ui/Button";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = {
      dataset:
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86658&format=file",
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
        <DataInput
          label="Data URL"
          value={this.state.dataset}
          onChange={this.handleChange("dataset")}
          margin="normal"
        />
        <Select
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
        </Select>
        <Spacer />
        <Button
          variant="raised"
          color="primary"
          onClick={() => this.props.run(this.state)}
          disabled={!this.canRun()}
        >
          Run
        </Button>
      </Container>
    );
  }

  canRun() {
    const dataValid = this.state.dataset !== "";
    const alignerValid = Object.keys(this.props.aligners).includes(
      this.state.aligner
    );
    return dataValid && alignerValid;
  }
}

const Container = styled.div`
  display: flex;
  align-items: center;
`;

const DataInput = styled(TextField)`
  flex: 1;
  max-width: 500px;
  margin-right: 20px !important;
`;

const Select = styled(TextField)`
  width: 200px;
  margin-right: 20px !important;
`;

const Spacer = styled.div`
  flex: 1;
`;
