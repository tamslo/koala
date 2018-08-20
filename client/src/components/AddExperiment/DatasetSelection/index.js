import React, { Component } from "react";
import styled from "styled-components";
import MenuItem from "@material-ui/core/MenuItem";
import Divider from "@material-ui/core/Divider";
import InfoIcon from "@material-ui/icons/Info";
import Select from "../../mui-wrappers/inputs/Select";
import Loading from "../../Loading";
import AddDataset from "./AddDataset";
import DatasetInfo from "./DatasetInfo";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = { addingDataset: false, datasetInfo: null };
  }

  handleChange = event => {
    if (event.target.value === "add") {
      this.startAddingDataset();
    } else {
      this.props.setDataset(event.target.value);
    }
  };

  render() {
    const { addingDataset, datasetInfo } = this.state;
    const { dataset } = this.props;
    // Make sure status indicators such as isLoading are not dispayed
    const datasets = Object.keys(this.props.datasets).reduce(
      (actualDatasets, datasetId) => {
        // Dummy test for uuid
        return datasetId.split("-").length === 5
          ? { ...actualDatasets, [datasetId]: this.props.datasets[datasetId] }
          : actualDatasets;
      },
      {}
    );

    const loadingText =
      "Adding data set, this can take a while for big files...";
    return (
      <div>
        {this.props.datasets.areLoading && <Loading content={loadingText} />}
        <AddDataset
          addDataset={this.addDataset.bind(this)}
          cancel={this.cancelAddingDataset.bind(this)}
          open={addingDataset}
        />
        <DatasetInfo
          open={datasetInfo !== null}
          close={this.closeDatasetInfo.bind(this)}
          dataset={datasetInfo || null}
        />
        <Select
          label="Data Set"
          value={dataset || ""}
          onChange={this.handleChange}
        >
          {Object.keys(datasets).map(datasetId => (
            <DatasetItem key={datasetId} value={datasetId}>
              <DatasetName>{datasets[datasetId].name}</DatasetName>
              <StyledInfoIcon
                onClick={this.showDatasetInfo(datasets[datasetId])}
              />
            </DatasetItem>
          ))}
          <Divider />
          <MenuItem value="add">Add data set</MenuItem>
        </Select>
      </div>
    );
  }

  startAddingDataset() {
    this.setState({ addingDataset: true });
  }

  cancelAddingDataset() {
    this.setState({ addingDataset: false });
  }

  addDataset(dataset) {
    this.setState({ addingDataset: false }, () =>
      this.props.addDataset(dataset)
    );
  }

  showDatasetInfo = dataset => event => {
    event.stopPropagation(); // Stop process of adding dataset to experiment
    this.setState({ datasetInfo: dataset });
  };

  closeDatasetInfo() {
    this.setState({ datasetInfo: null });
  }
}

const StyledInfoIcon = styled(InfoIcon)`
  height: 20px !important;
  width: 20px !important;
  color: #666666;
  margin-left: 6px;
  &:hover {
    color: #333333;
  }
`;

const DatasetItem = styled(MenuItem)`
  justify-content: space-between !important;
`;

const DatasetName = styled.div`
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
`;
