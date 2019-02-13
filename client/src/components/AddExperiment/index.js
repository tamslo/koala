import React, { Component } from "react";
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
  VARIANT_FILTERING,
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
      name: "New Experiment",
      reference: "",
      datasets: [],
      pipeline: {
        [ALIGNMENT]: "",
        [ALIGNMENT_FILTERING]: [],
        [VARIANT_CALLING]: "",
        [VARIANT_FILTERING]: ""
      }
    };
  }

  canRun() {
    return (
      this.state.name !== "" &&
      this.state.reference !== "" &&
      this.state.datasets.length > 0 &&
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
          selected={this.state.datasets}
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

        <Select
          label={handlerName(VARIANT_FILTERING)}
          value={this.state.pipeline[VARIANT_FILTERING]}
          onChange={this.handlePipelineChange(VARIANT_FILTERING)}
        >
          <MenuItem key="no-filter" value="">
            <em>None</em>
          </MenuItem>
          {this.props.services
            .filter(service => service.type === serviceTypes[VARIANT_FILTERING])
            .map(filter => (
              <MenuItem key={filter.id} value={filter.id}>
                {filter.name}
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
    this.setState({ datasets: [...this.state.datasets, dataset.id] }, () =>
      this.props.addDataset(dataset)
    );
  }

  setDataset(datasetIds) {
    this.setState({ datasets: datasetIds });
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
