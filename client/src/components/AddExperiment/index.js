import React, { Component } from "react";
import styled from "styled-components";
import TextField from "material-ui/TextField";
import MenuItem from "material-ui/Menu/MenuItem";
import Button from "material-ui/Button";
import Divider from "material-ui/Divider";
import AddDataset from "./AddDataset";

export default class extends Component {
  constructor(props) {
    super(props);
    const aligners = this.props.services.filter(
      service => service.type === "aligner"
    );
    this.state = {
      name: "New Experiment",
      datasets: this.props.datasets,
      dataset: "",
      aligner: "",
      aligners,
      addDataset: false
    };
  }

  handleChange = name => event => {
    if (name === "dataset" && event.target.value === "add") {
      this.startAddDataset();
    } else {
      this.setState({
        [name]: event.target.value
      });
    }
  };

  render() {
    return (
      <Container>
        <AddDataset
          addDataset={this.addDataset.bind(this)}
          cancel={this.cancelAddDataset.bind(this)}
          open={this.state.addDataset}
        />
        <FixedWidthTextField
          label="Name"
          value={this.state.name}
          onChange={this.handleChange("name")}
          margin="normal"
        />
        <FixedWidthTextField
          select
          label="Data Set"
          value={this.state.dataset || ""}
          onChange={this.handleChange("dataset")}
          margin="normal"
        >
          {Object.keys(this.state.datasets).map(datasetId => (
            <MenuItem key={datasetId} value={datasetId}>
              {datasetId}
            </MenuItem>
          ))}
          <Divider />
          <MenuItem value="add">Add data set</MenuItem>
        </FixedWidthTextField>
        <FixedWidthTextField
          select
          label="Aligner"
          value={this.state.aligner || ""}
          onChange={this.handleChange("aligner")}
          margin="normal"
        >
          {this.state.aligners.map(aligner => (
            <MenuItem key={aligner.id} value={aligner.id}>
              {aligner.name}
            </MenuItem>
          ))}
        </FixedWidthTextField>
        <VerticalSpacer />
        <Button
          color="primary"
          onClick={() =>
            this.props.addExperiment({
              name: this.state.name,
              dataset: this.props.dataset(this.state.dataset),
              alignment: this.state.aligner
            })
          }
          disabled={!this.canRun()}
          size="large"
        >
          Add
        </Button>
      </Container>
    );
  }

  startAddDataset() {
    this.setState({ addDataset: true });
  }

  cancelAddDataset() {
    this.setState({ addDataset: false });
  }

  addDataset(datasetId, dataset) {
    this.setState({
      addDataset: false,
      datasets: { ...this.state.datasets, [datasetId]: dataset },
      dataset: datasetId
    });
  }

  canRun() {
    const nameValid = this.state.name !== "";
    const dataValid = this.state.dataset !== "";
    const alignerValid = this.state.aligners
      .map(aligner => aligner.id)
      .includes(this.state.aligner);
    return nameValid && dataValid && alignerValid;
  }
}

const Container = styled.div`
  display: flex;
  align-items: center;
  flex-wrap: wrap;
  padding-left: 12px;
`;

const FixedWidthTextField = styled(TextField)`
  width: 200px;
  margin-right: 20px !important;
`;

const VerticalSpacer = styled.div`
  flex: 1;
`;
