import React, { Component } from "react";

const INTERVAL = 3000;

let refreshInterval;

export default WrappedComponent => {
  return class extends Component {
    render() {
      clearInterval(refreshInterval);

      if (this.props.context !== null) {
        this.updateExperiment();
      }

      return <WrappedComponent {...this.props} />;
    }

    updateExperiment() {
      const { context, updateExperiment } = this.props;
      const runningExperimentIndex = Object.values(
        context.experiments
      ).findIndex(experiment => experiment.running);
      const runningExperiment =
        runningExperimentIndex >= 0 &&
        Object.keys(context.experiments)[runningExperimentIndex];
      if (runningExperiment) {
        refreshInterval = setInterval(() => {
          updateExperiment(runningExperiment);
        }, INTERVAL);
      }
    }
  };
};
