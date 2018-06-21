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
      content: {}
    };
  }

  canAdd() {
    const contentLength =
      this.state.layout === constants.dataset.PAIRED ? 2 : 1;
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
    if (layout === constants.dataset.SINGLE) {
      content = content[constants.dataset.FORWARD]
        ? { [constants.dataset.FORWARD]: content[constants.dataset.FORWARD] }
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
      this.state.method === constants.dataset.FILE
        ? event.target.files[0]
        : { name: event.target.value };
    const removeKey =
      this.state.method === constants.dataset.FILE
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
