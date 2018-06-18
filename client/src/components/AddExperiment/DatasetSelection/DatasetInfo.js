import React, { Component } from "react";
import Dialog from "../../mui-wrappers/Dialog";
import DatasetInputs from "./DatasetInputs";

export default class extends Component {
  render() {
    const { dataset, open, close } = this.props;
    const actions = [
      {
        name: "Close",
        onClick: close,
        color: "primary"
      }
    ];

    if (!dataset) {
      return null;
    }

    return (
      <Dialog open={open} title="Data Set Info" actions={actions}>
        <DatasetInputs {...dataset} disabled={true} />
      </Dialog>
    );
  }
}
