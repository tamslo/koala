import React, { Component } from "react";
import { connect } from "react-redux";

import { fetchContext } from "../actions";
import {
  addExperiment,
  deleteExperiment,
  runExperiment,
  retryExperiment,
  refreshExperiment
} from "../actions/experiment";
import { addDataset } from "../actions/dataset";

export default WrappedComponent => {
  class StoreWrapper extends Component {
    componentWillMount() {
      this.props.fetchContext();
    }

    render() {
      return <WrappedComponent {...this.props} />;
    }
  }

  const mapStateToProps = state => {
    return {
      context: state.context,
      jobs: state.jobs
    };
  };

  const actions = {
    fetchContext,
    addExperiment,
    deleteExperiment,
    runExperiment,
    retryExperiment,
    refreshExperiment,
    addDataset
  };

  return connect(
    mapStateToProps,
    actions
  )(StoreWrapper);
};
