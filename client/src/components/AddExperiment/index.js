import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import MenuItem from "@material-ui/core/MenuItem";
import Button from "@material-ui/core/Button";
import TextField from "../mui-wrappers/inputs/Text";
import Select from "../mui-wrappers/inputs/Select";
import DatasetSelection from "./DatasetSelection";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      name: "New Experiment",
      reference: "",
      dataset: "",
      alignment: ""
    };
  }

  canRun() {
    return (
      this.state.name !== "" &&
      this.state.reference !== "" &&
      this.state.dataset !== "" &&
      this.state.alignment !== ""
    );
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  render() {
    return (
      <Container>
        <TextField
          label="Name"
          value={this.state.name}
          onChange={this.handleChange("name")}
        />

        <Select
          label="Reference Genome"
          value={this.state.reference || ""}
          onChange={this.handleChange("reference")}
        >
          {this.props.services
            .filter(service => service.type === "reference")
            .map(genome => (
              <MenuItem key={genome.id} value={genome.id}>
                {genome.name}
              </MenuItem>
            ))}
        </Select>

        <DatasetSelection
          datasets={this.props.datasets}
          dataset={this.state.dataset}
          addDataset={this.addDataset.bind(this)}
          setDataset={this.setDataset.bind(this)}
        />

        <Select
          label="Aligner"
          value={this.state.alignment || ""}
          onChange={this.handleChange("alignment")}
        >
          {this.props.services
            .filter(service => service.type === "aligner")
            .map(aligner => (
              <MenuItem key={aligner.id} value={aligner.id}>
                {aligner.name}
              </MenuItem>
            ))}
        </Select>

        <VerticalSpacer />
        <Button
          color="primary"
          onClick={this.addExperiment.bind(this)}
          disabled={!this.canRun()}
          size="large"
        >
          Add
        </Button>
      </Container>
    );
  }

  addExperiment() {
    const experiment = this.state;
    this.setState(this.initialState(), () =>
      this.props.addExperiment(experiment)
    );
  }

  addDataset(dataset) {
    this.setState({ dataset: dataset.id }, () =>
      this.props.addDataset(dataset)
    );
  }

  setDataset(dataset) {
    this.setState({ dataset });
  }
}

const Container = styled.div`
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  padding-left: 12px;
`;

const VerticalSpacer = styled.div`
  flex: 1;
`;
