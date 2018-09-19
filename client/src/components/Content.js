import React, { Component } from "react";
import styled from "styled-components";
import Loading from "./Loading";
import AddExperiment from "./AddExperiment";
import Experiments from "./Experiments";

export default class extends Component {
  render() {
    const {
      context,
      addExperiment,
      deleteExperiment,
      retryExperiment,
      addDataset
    } = this.props;

    if (context === null) {
      return <Loading content={"Setting everything up..."} />;
    }

    if (context.error) {
      return (
        <Loading
          content={context.error.message}
          error={context.error}
          retry={this.props.fetchContext}
        />
      );
    }

    return (
      <Container className="content">
        <AddExperiment
          services={context.services}
          addExperiment={addExperiment}
          datasets={context.datasets}
          addDataset={addDataset}
        />
        <Spacer />
        {Object.keys(context.experiments).length > 0 && (
          <Experiments
            experiments={context.experiments}
            services={context.services}
            datasets={context.datasets}
            deleteExperiment={deleteExperiment}
            retryExperiment={retryExperiment}
          />
        )}
      </Container>
    );
  }
}

const Container = styled.div`
  padding: 32px;
`;

const Spacer = styled.div`
  height: 32px;
`;
