import React, { Component } from "react";
import uuid from "uuid/v4";
import styled from "styled-components";
import TextField from "@material-ui/core/TextField";
import Dialog from "../mui-wrappers/Dialog";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return { id: uuid(), name: "", url: "" };
  }

  canAdd() {
    return this.state.name !== "" && this.state.url !== "";
  }

  handleChange = name => event => {
    this.setState({
      [name]: event.target.value
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
        onClick: () => this.props.addDataset(this.state),
        color: "primary",
        disabled: !this.canAdd()
      }
    ];

    return (
      <Dialog open={this.props.open} title="Add Data Set" actions={actions}>
        <StyledTextField
          label="Name"
          value={this.state.name}
          onChange={this.handleChange("name")}
          margin="normal"
        />
        <StyledTextField
          label="Data URL"
          value={this.state.url}
          onChange={this.handleChange("url")}
          margin="normal"
        />
      </Dialog>
    );
  }

  addDataset() {
    this.setState(this.initialState(), () => this.props.addDataset(this.state));
  }
}

const StyledTextField = styled(TextField)`
  flex-grow: 1;
  min-width: 200px;
  margin-right: 20px !important;
`;
