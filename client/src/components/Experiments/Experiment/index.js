import React, { Component } from "react";
import styled from "styled-components";
import IconButton from "@material-ui/core/IconButton";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import Log from "./Log";

export default class extends Component {
  render() {
    const dataset = this.props.datasets[this.props.dataset];
    return (
      <div>
        <Entry>
          {`Data Set: ${dataset.name}`}
          {this.renderDownloadDatasetButton(dataset)}
        </Entry>
        <Entry>
          {`Aligner: ${
            this.props.services.find(
              service => service.id === this.props.pipeline.alignment.id
            ).name
          }`}
          {this.renderDownloadButton(this.props.pipeline["alignment"].file)}
        </Entry>
        <Log
          pipeline={this.props.pipeline}
          error={this.props.error}
          created={this.props.created}
          done={this.props.done}
          interrupted={this.props.interrupted}
        />
      </div>
    );
  }

  renderDownloadDatasetButton(dataset) {
    const { data } = dataset;
    const firstFilePath = data[Object.keys(data)[0]].path;
    let path;
    if (Object.keys(data).length === 1) {
      path = firstFilePath;
    } else {
      const separator = "/";
      path = firstFilePath
        .split(separator)
        .slice(0, -1)
        .join(separator);
    }
    return this.renderDownloadButton(path);
  }

  renderDownloadButton(path) {
    return path ? (
      <StyledIconButton
        aria-label={"Download"}
        href={this.props.SERVER_URL + "/export?path=" + path}
      >
        <DownloadIcon />
      </StyledIconButton>
    ) : null;
  }
}

const Entry = styled.div`
  margin-bottom: 12px;
`;

const StyledIconButton = styled(IconButton)`
  margin-left: 12px !important;
`;
