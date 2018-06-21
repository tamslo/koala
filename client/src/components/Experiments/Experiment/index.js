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
          {`Data Set: ${this.props.datasets[this.props.dataset].name}`}
          {this.renderDownloadButton("dataset")}
        </Entry>
        <Entry>
          {`Aligner: ${
            this.props.services.find(
              service => service.id === this.props.alignment
            ).name
          }`}
          {this.renderDownloadButton("alignment")}
        </Entry>
        <Log log={this.props.log} error={this.props.error} />
      </div>
    );
  }

  renderDownloadButton(key) {
    const path = this.props.files[key];
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
