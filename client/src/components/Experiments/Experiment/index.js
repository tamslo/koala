import React, { Component } from "react";
import styled from "styled-components";
import IconButton from "@material-ui/core/IconButton";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import LinkIcon from "@material-ui/icons/Link";
import Log from "./Log";
import { displayNames } from "../../experimentConstants";
import { handlerName } from "../../experimentUtils";

export default class extends Component {
  render() {
    const dataset = this.props.datasets[this.props.dataset];
    const reference = this.props.references.find(
      reference => reference.id === this.props.reference
    );
    return (
      <div>
        <Entry>
          {`${displayNames.reference}: ${reference.name}`}
          <StyledIconButton aria-label={"Source"} href={reference.source}>
            <LinkIcon />
          </StyledIconButton>
        </Entry>
        <Entry>
          {`${displayNames.dataset}: ${dataset.name}`}
          {this.renderDownloadDatasetButton(dataset)}
        </Entry>
        {Object.keys(this.props.pipeline).map(
          this.renderPipelineEntry.bind(this)
        )}
        <Log
          pipeline={this.props.pipeline}
          created={this.props.created}
          status={this.props.status}
          services={this.props.services}
        />
      </div>
    );
  }

  renderPipelineEntry(actionName) {
    const service = this.props.services.find(
      service => service.id === this.props.pipeline[actionName].id
    );
    return (
      <Entry key={`pipeline-entry-${actionName}`}>
        {`${handlerName(actionName)}: ${service.name}`}
        {this.renderDownloadButton(this.props.pipeline[actionName].directory)}
      </Entry>
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
  line-height: 32px;
`;

const StyledIconButton = styled(IconButton)`
  margin-left: 12px !important;
  width: 32px !important;
  height: 32px !important;
  .MuiSvgIcon-root-144 {
    font-size: 18px !important;
  }
`;
