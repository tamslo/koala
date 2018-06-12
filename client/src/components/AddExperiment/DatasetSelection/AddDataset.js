import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import Button from "@material-ui/core/Button";
import MenuItem from "@material-ui/core/MenuItem";
import Dialog from "../../mui-wrappers/Dialog";
import TextField from "../../mui-wrappers/inputs/Text";
import NumberField from "../../mui-wrappers/inputs/Number";
import Select from "../../mui-wrappers/inputs/Select";
import Checkbox from "../../mui-wrappers/inputs/Checkbox";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      name: "",
      method: "file",
      content: "",
      pairedEnd: true,
      readLength: 200
    };
  }

  canAdd() {
    return (
      this.state.name !== "" &&
      this.state.content !== "" &&
      this.state.readLength !== ""
    );
  }

  handleChange = name => event => {
    let value;
    if (name === "pairedEnd") {
      value = event.target.checked;
    } else if (name === "content") {
      value = event.target.files[0];
    } else {
      value = event.target.value;
    }
    this.setState({
      [name]: value
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
          />
          <Select
            label="Method"
            value={this.state.method}
            onChange={this.handleChange("method")}
          >
            <MenuItem value="file">File</MenuItem>
            <MenuItem value="url">URL</MenuItem>
          </Select>
          {this.renderFileSelection()}
          <Checkbox
            label="Paired end"
            onChange={this.handleChange("pairedEnd")}
            checked={this.state.pairedEnd}
          />
          <NumberField
            label="Read length"
            onChange={this.handleChange("readLength")}
            value={this.state.readLength}
          />
        </Container>
      </Dialog>
    );
  }

  renderFileSelection() {
    return this.state.method === "file" ? (
      <StyledButton
        variant="outlined"
        onClick={() => {
          this.refs.file.click();
        }}
      >
        Upload file
        <input
          ref="file"
          type="file"
          style={{ display: "none" }}
          onChange={this.handleChange("content")}
        />
      </StyledButton>
    ) : (
      <TextField
        label="Data URL"
        value={this.state.url}
        onChange={this.handleChange("content")}
      />
    );
  }

  addDataset() {
    const dataset = this.state;
    this.setState(this.initialState(), () => this.props.addDataset(dataset));
  }
}

const Container = styled.div`
  display: flex;
  flex-wrap: wrap;
`;

const StyledButton = styled(Button)`
  margin-right: 20px !important;
`;
