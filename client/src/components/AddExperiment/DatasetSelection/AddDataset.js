import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import Button from "@material-ui/core/Button";
import MenuItem from "@material-ui/core/MenuItem";
import Dialog from "../../mui-wrappers/Dialog";
import TextField from "../../mui-wrappers/inputs/Text";
import NumberField from "../../mui-wrappers/inputs/Number";
import Select from "../../mui-wrappers/inputs/Select";

const FORWARD = "forward";
const REVERSE = "reverse";
const PAIRED = "paired";
const SINGLE = "single";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      name: "New Data Set",
      layout: PAIRED,
      readLength: 200,
      method: "file",
      content: {}
    };
  }

  canAdd() {
    const contentLength = this.state.layout === PAIRED ? 2 : 1;
    return (
      this.state.name !== "" &&
      Object.keys(this.state.content).length === contentLength &&
      this.state.readLength !== ""
    );
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  changeLayout = event => {
    const layout = event.target.value;
    let content = this.state.content;

    // Remove reverse file from content if layout changes to single end
    if (layout === SINGLE) {
      content = content[FORWARD]
        ? { [FORWARD]: content[FORWARD] }
        : this.initialState().content;
    }

    this.setState({
      layout,
      content
    });
  };

  changeMethod = event => {
    this.setState({
      method: event.target.value,
      content: this.initialState().content
    });
  };

  changeContent = key => event => {
    const value =
      this.state.method === "file" ? event.target.files[0] : event.target.value;
    const removeKey =
      this.state.method === "file"
        ? event.target.files.length === 0
        : event.target.value === "";
    const contentWithoutKey = Object.keys(this.state.content).reduce(
      (reducedContent, fileKey) =>
        fileKey === key
          ? reducedContent
          : { ...reducedContent, [fileKey]: this.state.content[fileKey] },
      {}
    );
    const contentWithNewValue = {
      ...this.state.content,
      [key]: value
    };
    this.setState({
      content: removeKey ? contentWithoutKey : contentWithNewValue
    });
  };

  render() {
    const actions = [
      {
        name: "Cancel",
        onClick: this.props.cancel
      },
      {
        name: "Add",
        onClick: this.addDataset.bind(this),
        color: "primary",
        disabled: !this.canAdd()
      }
    ];

    return (
      <Dialog open={this.props.open} title="Add Data Set" actions={actions}>
        <Container>
          <TextField
            label="Name"
            value={this.state.name}
            onChange={this.handleChange("name")}
            width={390}
          />
          <Row>
            <NumberField
              label="Read length"
              onChange={this.handleChange("readLength")}
              value={this.state.readLength}
              width={100}
            />
            <Select
              label="Method"
              value={this.state.method}
              onChange={this.changeMethod}
              width={100}
            >
              <MenuItem value="file">File</MenuItem>
              <MenuItem value="url">URL</MenuItem>
            </Select>
            <Select
              label="Layout"
              value={this.state.layout}
              onChange={this.changeLayout}
              width={150}
            >
              <MenuItem value={PAIRED}>Paired end</MenuItem>
              <MenuItem value={SINGLE}>Single end</MenuItem>
            </Select>
          </Row>
          {this.renderDataSelection()}
        </Container>
      </Dialog>
    );
  }

  renderDataSelection() {
    let selections = [this.renderSingleDataSelection(FORWARD)];
    if (this.state.layout === PAIRED) {
      selections = [...selections, this.renderSingleDataSelection(REVERSE)];
    }
    return selections;
  }

  renderSingleDataSelection(key) {
    return this.state.method === "file" ? (
      this.renderFileUpload(key)
    ) : (
      <TextField
        key={key}
        label={this.label("Data URL", key)}
        value={this.state.url}
        onChange={this.changeContent(key)}
        width={390}
      />
    );
  }

  renderFileUpload(key) {
    const fileName =
      this.state.content[key] &&
      typeof this.state.content[key] === "object" &&
      this.state.content[key].name;
    const label = this.label(fileName || "Select file", key);
    return (
      <div key={key}>
        <StyledButton
          variant="outlined"
          onClick={() => {
            this.refs[key].click();
          }}
        >
          {label}
          <input
            ref={key}
            type="file"
            style={{ display: "none" }}
            onChange={this.changeContent(key)}
          />
        </StyledButton>
      </div>
    );
  }

  addDataset() {
    const dataset = this.state;
    this.setState(this.initialState(), () => this.props.addDataset(dataset));
  }

  label(text, key) {
    if (this.state.layout === PAIRED) {
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
