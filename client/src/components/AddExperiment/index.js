import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import MenuItem from "@material-ui/core/MenuItem";
import Button from "@material-ui/core/Button";
import TextField from "../mui-wrappers/inputs/Text";
import Select from "../mui-wrappers/inputs/Select";
import MultipleSelect from "../mui-wrappers/inputs/MultipleSelect";
import DatasetSelection from "./DatasetSelection";
import {
  ALIGNMENT,
  ALIGNMENT_FILTERING,
  VARIANT_CALLING,
  serviceTypes,
  displayNames
} from "../experimentConstants";
import { handlerName } from "../experimentUtils";

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
      pipeline: {
        [ALIGNMENT]: "",
        [ALIGNMENT_FILTERING]: [],
        [VARIANT_CALLING]: ""
      }
    };
  }

  canRun() {
    return (
      this.state.name !== "" &&
      this.state.reference !== "" &&
      this.state.dataset !== "" &&
      this.state.pipeline[ALIGNMENT] !== ""
    );
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  handlePipelineChange = name => event => {
    this.setState({
      pipeline: { ...this.state.pipeline, [name]: event.target.value }
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
          label={displayNames.reference}
          value={this.state.reference || ""}
          onChange={this.handleChange("reference")}
        >
          {this.props.references.map(genome => (
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
          label={handlerName(ALIGNMENT)}
          value={this.state.pipeline[ALIGNMENT]}
          onChange={this.handlePipelineChange(ALIGNMENT)}
        >
          {this.props.services
            .filter(service => service.type === serviceTypes[ALIGNMENT])
            .map(aligner => (
              <MenuItem key={aligner.id} value={aligner.id}>
                {aligner.name}
              </MenuItem>
            ))}
        </Select>

        <MultipleSelect
          id="alignment-filters"
          label={handlerName(ALIGNMENT_FILTERING)}
          selected={this.state.pipeline[ALIGNMENT_FILTERING]}
          onChange={this.handlePipelineChange(ALIGNMENT_FILTERING)}
          items={this.props.services.filter(
            service => service.type === serviceTypes[ALIGNMENT_FILTERING]
          )}
        />

        <Select
          label={handlerName(VARIANT_CALLING)}
          value={this.state.pipeline[VARIANT_CALLING]}
          onChange={this.handlePipelineChange(VARIANT_CALLING)}
        >
          <MenuItem key="no-caller" value="">
            <em>None</em>
          </MenuItem>
          {this.props.services
            .filter(service => service.type === serviceTypes[VARIANT_CALLING])
            .map(caller => (
              <MenuItem key={caller.id} value={caller.id}>
                {caller.name}
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

  setDataset(datasetId) {
    this.setState({ dataset: datasetId });
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
