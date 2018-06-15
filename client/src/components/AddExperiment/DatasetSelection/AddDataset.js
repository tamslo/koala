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
      name: "New Data Set",
      pairedEnd: true,
      readLength: 200,
      method: "file",
      content: {}
    };
  }

  canAdd() {
    // TODO check that URL is not empty and valid
    return (
      this.state.name !== "" &&
      Object.keys(this.state.content).length !== "" &&
      this.state.readLength !== ""
    );
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
    });
  };

  changePairedEnd = event => {
    // TODO remove file from content if true changes to false
    this.setState({
      pairedEnd: event.target.checked
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
    this.setState({
      content: {
        ...this.state.content,
        [key]: value
      }
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
            width={400}
          />
          <Row>
            <Select
              label="Method"
              value={this.state.method}
              onChange={this.changeMethod}
              width={50}
            >
              <MenuItem value="file">File</MenuItem>
              <MenuItem value="url">URL</MenuItem>
            </Select>
            {this.renderFileSelection()}
          </Row>
          <div>
            <NumberField
              label="Read length"
              onChange={this.handleChange("readLength")}
              value={this.state.readLength}
              width={100}
            />
            <Checkbox
              label="Paired end"
              onChange={this.changePairedEnd}
              checked={this.state.pairedEnd}
            />
          </div>
        </Container>
      </Dialog>
    );
  }

  renderFileSelection() {
    return this.state.method === "file" ? (
      this.renderFileUpload()
    ) : (
      <TextField
        label="Data URL"
        value={this.state.url}
        onChange={this.changeContent("forward")}
        width={330}
      />
    );
  }

  renderFileUpload() {
    const fileName =
      this.state.content.forward &&
      typeof this.state.content.forward === "object" &&
      this.state.content.forward.name;
    return (
      <div>
        <Button
          variant="outlined"
          onClick={() => {
            this.refs.file.click();
          }}
        >
          {fileName || "Select file"}
          <input
            ref="file"
            type="file"
            style={{ display: "none" }}
            onChange={this.changeContent("forward")}
          />
        </Button>
      </div>
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
  flex-direction: column;
`;

const Row = styled.div`
  display: flex;
  align-items: baseline;
`;
