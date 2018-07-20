import React, { Component } from "react";
import uuid from "uuid/v4";
import Dialog from "../../mui-wrappers/Dialog";
import DatasetInputs from "./DatasetInputs";
import constants from "../../../constants.json";

export default class extends Component {
  constructor(props) {
    super(props);
    this.state = this.initialState();
  }

  initialState() {
    return {
      id: uuid(),
      name: "New Data Set",
      layout: constants.dataset.PAIRED,
      readLength: 200,
      method: constants.dataset.FILE,
      data: {}
    };
  }

  canAdd() {
    const dataLength = this.state.layout === constants.dataset.PAIRED ? 2 : 1;
    return (
      this.state.name !== "" &&
      Object.keys(this.state.data).length === dataLength &&
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
    let data = this.state.data;

    // Remove reverse file from data if layout changes to single end
    if (layout === constants.dataset.SINGLE) {
      data = data[constants.dataset.FORWARD]
        ? { [constants.dataset.FORWARD]: data[constants.dataset.FORWARD] }
        : this.initialState().data;
    }

    this.setState({
      layout,
      data
    });
  };

  changeMethod = event => {
    this.setState({
      method: event.target.value,
      data: this.initialState().data
    });
  };

  changeContent = key => event => {
    const value =
      this.state.method === constants.dataset.FILE
        ? event.target.files[0]
        : { name: event.target.value };
    const removeKey =
      this.state.method === constants.dataset.FILE
        ? event.target.files.length === 0
        : event.target.value === "";
    const dataWithoutKey = Object.keys(this.state.data).reduce(
      (reducedContent, fileKey) =>
        fileKey === key
          ? reducedContent
          : { ...reducedContent, [fileKey]: this.state.data[fileKey] },
      {}
    );
    const dataWithNewValue = {
      ...this.state.data,
      [key]: value
    };
    this.setState({
      data: removeKey ? dataWithoutKey : dataWithNewValue
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
        <DatasetInputs
          {...this.state}
          handleChange={this.handleChange.bind(this)}
          changeLayout={this.changeLayout.bind(this)}
          changeMethod={this.changeMethod.bind(this)}
          changeContent={this.changeContent.bind(this)}
        />
      </Dialog>
    );
  }

  addDataset() {
    const dataset = this.state;
    this.setState(this.initialState(), () => this.props.addDataset(dataset));
  }
}
