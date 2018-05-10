import React, { Component } from "react";
import { connect } from "react-redux";
import styled from "styled-components";
import { fetchContext, run } from "../actions";
import Loading from "./Loading";
import Inputs from "./Inputs";
import Experiments from "./Experiments";

class Content extends Component {
  componentWillMount() {
    this.props.fetchContext();
  }

  render() {
    const { context } = this.props;
    if (context === null) {
      return <Loading content={"Setting everything up..."} />;
    }

    if (context.isError) {
      return (
        <Loading
          content={context.error.message}
          error={context.isError}
          retry={this.props.fetchContext}
        />
      );
    }

    return (
      <Container className="content">
        <Inputs aligners={this.props.context.aligners} run={this.props.run} />
        {Object.keys(this.props.context.experiments).length > 0 && (
          <Experiments experiments={this.props.context.experiments} />
        )}
      </Container>
    );
  }
}

const Container = styled.div`
  padding: 10px;
`;

const mapStateToProps = state => {
  return {
    context: state
  };
};

const mapDispatchToProps = dispatch => {
  return {
    fetchContext: () => {
      dispatch(fetchContext());
    },
    run: params => {
      dispatch(run(params));
    }
  };
};

export default connect(mapStateToProps, mapDispatchToProps)(Content);
