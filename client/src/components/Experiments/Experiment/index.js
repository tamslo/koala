import React, { Component } from "react";
import styled from "styled-components";
import IconButton from "@material-ui/core/IconButton";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import LinkIcon from "@material-ui/icons/Link";
import Log from "./Log";

export default class extends Component {
  render() {
    const dataset = this.props.datasets[this.props.dataset];
    const reference = this.getWithId(
      this.props.references,
      this.props.reference
    );
    const aligner = this.getWithId(
      this.props.services,
      this.props.pipeline.alignment.id
    );
    return (
      <div>
        <Entry>
          {`Reference Genome: ${reference.name}`}
          <StyledIconButton aria-label={"Source"} href={reference.source}>
            <LinkIcon />
          </StyledIconButton>
        </Entry>
        <Entry>
          {`Data Set: ${dataset.name}`}
          {this.renderDownloadDatasetButton(dataset)}
        </Entry>
        <Entry>
          {`Aligner: ${aligner.name}`}
          {this.renderDownloadButton(this.props.pipeline["alignment"].file)}
        </Entry>
        <Log
          pipeline={this.props.pipeline}
          created={this.props.created}
          status={this.props.status}
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

  getWithId(collection, id) {
    return collection.find(item => item.id === id);
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
