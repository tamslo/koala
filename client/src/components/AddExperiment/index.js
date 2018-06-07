import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import TextField from "material-ui/TextField";
import MenuItem from "material-ui/Menu/MenuItem";
import Button from "material-ui/Button";
import Divider from "material-ui/Divider";
import InfoIcon from "@material-ui/icons/Info";
import AddDataset from "./AddDataset";
import DatasetInfo from "./DatasetInfo";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      addDataset: false,
      datasetInfo: null,
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
    const { dataset, addDataset, datasetInfo, name, aligner } = this.state;
    const { datasets, services } = this.props;

    return (
      <Container>
        <AddDataset
          addDataset={this.addDataset.bind(this)}
          cancel={this.cancelAddDataset.bind(this)}
          open={addDataset}
        />
        <DatasetInfo
          open={datasetInfo !== null}
          close={this.closeDatasetInfo.bind(this)}
          dataset={(datasetInfo && datasets[dataset]) || {}}
        />
        <FixedWidthTextField
          label="Name"
          value={name}
          onChange={this.handleChange("name")}
          margin="normal"
        />
        {this.renderDatasetSelection()}
        <FixedWidthTextField
          select
          label="Aligner"
          value={aligner || ""}
          onChange={this.handleChange("aligner")}
          margin="normal"
        >
          {services
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
          onClick={this.addExperiment.bind(this)}
          disabled={!this.canRun()}
          size="large"
        >
          Add
        </Button>
      </Container>
    );
  }

  renderDatasetSelection() {
    const { dataset } = this.state;
    const { datasets } = this.props;
    return (
      <FixedWidthTextField
        select
        label="Data Set"
        value={dataset || ""}
        onChange={this.handleChange("dataset")}
        margin="normal"
      >
        {Object.keys(datasets).map(datasetId => (
          <DatasetItem key={datasetId} value={datasetId}>
            {datasets[datasetId].name}
            <StyledInfoIcon
              onClick={event => {
                this.showDatasetInfo(datasets[datasetId], event);
              }}
            />
          </DatasetItem>
        ))}
        <Divider />
        <MenuItem value="add">Add data set</MenuItem>
      </FixedWidthTextField>
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

  startAddDataset() {
    this.setState({ addDataset: true });
  }

  cancelAddDataset() {
    this.setState({ addDataset: false });
  }

  addDataset(dataset) {
    this.setState({ addDataset: false, dataset: dataset.id }, () =>
      this.props.addDataset(dataset)
    );
  }

  showDatasetInfo(dataset, event) {
    console.log(event.target);
    console.log(this.state);
    this.setState({ datasetInfo: dataset, dataset: "" });
  }

  closeDatasetInfo() {
    this.setState({ datasetInfo: null });
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
  justify-content: space-between !important;
`;

const VerticalSpacer = styled.div`
  flex: 1;
`;

const StyledInfoIcon = styled(InfoIcon)`
  height: 20px !important;
  width: 20px !important;
  color: #666666;
`;

const DatasetItem = styled(MenuItem)`
  justify-content: space-between !important;
`;
