import React, { Component } from "react";

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
          const experimentWaiting =
            !experiment.done && !experiment.error && !experiment.interrupted;
          return experiment.running || experimentWaiting;
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
