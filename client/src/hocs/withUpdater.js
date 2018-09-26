import React, { Component } from "react";
import constants from "../constants.json";

const INTERVAL = 3000;

let refreshInterval;

export default WrappedComponent => {
  return class extends Component {
    render() {
      clearInterval(refreshInterval);

      if (this.props.context !== null) {
        this.updateRunningExperiment();
      }

      return <WrappedComponent {...this.props} />;
    }

    updateRunningExperiment() {
      const { context, updateRunningExperiment } = this.props;
      const shouldUpdate = Object.values(context.experiments).some(
        experiment => {
          return (
            experiment.status === constants.experiment.RUNNING ||
            experiment.status === constants.experiment.WAITING
          );
        }
      );
      if (shouldUpdate) {
        refreshInterval = setInterval(() => {
          updateRunningExperiment();
        }, INTERVAL);
      }
    }
  };
};
