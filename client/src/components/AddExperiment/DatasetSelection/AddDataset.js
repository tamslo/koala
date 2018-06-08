import React, { Component } from "react";
import uuid from "uuid/v4";
import Dialog from "../../mui-wrappers/Dialog";
import TextField from "../../mui-wrappers/TextField";

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
        <TextField
          label="Name"
          value={this.state.name}
          onChange={this.handleChange("name")}
        />
        <TextField
          label="Data URL"
          value={this.state.url}
          onChange={this.handleChange("url")}
        />
      </Dialog>
    );
  }

  addDataset() {
    this.setState(this.initialState(), () => this.props.addDataset(this.state));
  }
}
