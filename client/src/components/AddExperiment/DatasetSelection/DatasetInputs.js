import React, { Component } from "react";
import styled from "styled-components";
import Button from "@material-ui/core/Button";
import MenuItem from "@material-ui/core/MenuItem";
import TextField from "../../mui-wrappers/inputs/Text";
import NumberField from "../../mui-wrappers/inputs/Number";
import Select from "../../mui-wrappers/inputs/Select";
import constants from "../../../constants.json";

export default class extends Component {
  render() {
    return (
      <Container>
        <TextField
          label="Name"
          value={this.props.name}
          onChange={this.props.handleChange && this.props.handleChange("name")}
          width={390}
          disabled={this.props.disabled}
        />
        <Row>
          <NumberField
            label="Read length"
            onChange={
              this.props.handleChange && this.props.handleChange("readLength")
            }
            value={this.props.readLength}
            width={100}
            disabled={this.props.disabled}
          />
          <Select
            label="Method"
            value={this.props.method}
            onChange={this.props.changeMethod}
            width={100}
            disabled={this.props.disabled}
          >
            <MenuItem value={constants.dataset.FILE}>File</MenuItem>
            <MenuItem value={constants.dataset.URL}>URL</MenuItem>
          </Select>
          <Select
            label="Layout"
            value={this.props.layout}
            onChange={this.props.changeLayout}
            width={150}
            disabled={this.props.disabled}
          >
            <MenuItem value={constants.dataset.PAIRED}>Paired end</MenuItem>
            <MenuItem value={constants.dataset.SINGLE}>Single end</MenuItem>
          </Select>
        </Row>
        {this.renderDataSelection()}
      </Container>
    );
  }

  renderDataSelection() {
    let selections = [
      this.renderSingleDataSelection(constants.dataset.FORWARD)
    ];
    if (this.props.layout === constants.dataset.PAIRED) {
      selections = [
        ...selections,
        this.renderSingleDataSelection(constants.dataset.REVERSE)
      ];
    }
    return selections;
  }

  renderSingleDataSelection(key) {
    const value =
      this.props.data[key] &&
      typeof this.props.data[key] === "object" &&
      this.props.data[key].name;
    return this.props.method === constants.dataset.FILE
      ? this.renderFileUpload(key, value)
      : this.renderUrlInput(key, value);
  }

  renderFileUpload(key, name) {
    const label = name || this.label("Select file", key);
    return (
      <div key={key}>
        <StyledButton
          variant="outlined"
          onClick={() => {
            this.refs[key].click();
          }}
          disabled={this.props.disabled}
        >
          {label}
          <input
            ref={key}
            type="file"
            style={{ display: "none" }}
            onChange={this.props.changeContent && this.props.changeContent(key)}
          />
        </StyledButton>
      </div>
    );
  }

  renderUrlInput(key, url = "") {
    return (
      <TextField
        key={key}
        label={this.label("Data URL", key)}
        value={url}
        onChange={this.props.changeContent && this.props.changeContent(key)}
        width={390}
        disabled={this.props.disabled}
      />
    );
  }

  label(text, key) {
    if (this.props.layout === constants.dataset.PAIRED) {
      text += ` (${key})`;
    }
    return text;
  }
}

const Container = styled.div`
  display: flex;
  flex-wrap: wrap;
  flex-direction: column;
`;

const Row = styled.div`
  display: flex;
  flex-direction: row;
`;

const StyledButton = styled(Button)`
  margin-top: 20px !important;
`;
