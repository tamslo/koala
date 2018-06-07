import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import TextField from "material-ui/TextField";
import MenuItem from "material-ui/Menu/MenuItem";
import Button from "material-ui/Button";
import Divider from "material-ui/Divider";
import AddDataset from "./AddDataset";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = {
      id: uuid(),
      addDataset: false,
      name: "",
      aligner: "",
      dataset: ""
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
          {Object.keys(this.props.datasets).map(datasetId => (
            <MenuItem key={datasetId} value={datasetId}>
              {this.props.datasets[datasetId].name}
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
          {this.props.services
            .filter(service => service.type === "aligner")
            .map(aligner => (
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
              id: this.state.id,
              name: this.state.name,
              dataset: this.props.datasets[this.state.dataset].id,
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

  addDataset(dataset) {
    this.setState(
      { addDataset: false, dataset: dataset.id },
      this.props.addDataset(dataset)
    );
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

const FixedWidthTextField = styled(TextField)`
  width: 200px;
  margin-right: 20px !important;
`;

const VerticalSpacer = styled.div`
  flex: 1;
`;
