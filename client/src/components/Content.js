import React, { Component } from "react";
import { connect } from "react-redux";
import styled from "styled-components";
import { fetchContext, run } from "../actions";
import LoadingDialog from "./LoadingDialog";
import Inputs from "./Inputs";

class Content extends Component {
  componentWillMount() {
    this.props.fetchContext();
  }

  render() {
    const { context } = this.props;
    if (context === null) {
      return <LoadingDialog content={"Setting everything up..."} />;
    }

    if (context.isError) {
      return (
        <LoadingDialog
          content={context.error.message}
          error={context.isError}
          retry={this.props.fetchContext}
        />
      );
    }

    return (
      <Container className="content">
        <Inputs aligners={this.props.context.aligners} run={this.props.run} />
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
