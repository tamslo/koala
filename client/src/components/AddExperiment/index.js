import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import MenuItem from "@material-ui/core/MenuItem";
import Button from "@material-ui/core/Button";
import TextField from "../mui-wrappers/TextField";
import DatasetSelection from "./DatasetSelection";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      name: "",
      aligner: "",
      dataset: ""
    };
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

        <DatasetSelection
          datasets={this.props.datasets}
          dataset={this.state.dataset}
          addDataset={this.addDataset.bind(this)}
          setDataset={this.setDataset.bind(this)}
        />

        <TextField
          select={true}
          label="Aligner"
          value={this.state.aligner || ""}
          onChange={this.handleChange("aligner")}
        >
          {this.props.services
            .filter(service => service.type === "aligner")
            .map(aligner => (
              <MenuItem key={aligner.id} value={aligner.id}>
                {aligner.name}
              </MenuItem>
            ))}
        </TextField>

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
    this.setState(
      this.initialState(),
      this.props.addExperiment({
        id: this.state.id,
        name: this.state.name,
        dataset: this.props.datasets[this.state.dataset].id,
        alignment: this.state.aligner
      })
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

  canRun() {
    return (
      this.state.name !== "" &&
      this.state.dataset !== "" &&
      this.state.aligner !== ""
    );
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
