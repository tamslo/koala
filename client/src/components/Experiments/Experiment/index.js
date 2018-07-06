import React, { Component } from "react";
import styled from "styled-components";
import IconButton from "@material-ui/core/IconButton";
import DownloadIcon from "@material-ui/icons/CloudDownload";
import Log from "./Log";

export default class extends Component {
  render() {
    return (
      <div>
        <Entry>
          {`Data Set: ${
            this.props.datasets[this.props.pipeline.dataset.id].name
          }`}
          {this.renderDownloadButton("dataset")}
        </Entry>
        <Entry>
          {`Aligner: ${
            this.props.services.find(
              service => service.id === this.props.pipeline.alignment.id
            ).name
          }`}
          {this.renderDownloadButton("alignment")}
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

  renderDownloadButton(action) {
    const path = this.props.pipeline[action].file;
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
