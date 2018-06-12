import React, { Component } from "react";
import uuid from "uuid/v4";
import Dialog from "../../mui-wrappers/Dialog";
import TextField from "../../mui-wrappers/TextField";
import Checkbox from "../../mui-wrappers/Checkbox";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return { id: uuid(), name: "", url: "", pairedEnd: true };
  }

  canAdd() {
    return this.state.name !== "" && this.state.url !== "";
  }

  handleChange = name => event => {
    let value;
    if (name === "pairedEnd") {
      value = event.target.checked;
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
        <Checkbox
          label="Paired End"
          onChange={this.handleChange("pairedEnd")}
          checked={this.state.pairedEnd}
        />
      </Dialog>
    );
  }

  addDataset() {
    const dataset = this.state;
    this.setState(this.initialState(), () => this.props.addDataset(dataset));
  }
}
