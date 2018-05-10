import React from "react";
import styled from "styled-components";
import Card, { CardContent } from "material-ui/Card";
import Typography from "material-ui/Typography";

export default props => {
  return (
    <StyledCard>
      <CardContent>
        <Typography variant="title">{props.title}</Typography>
        <Spacer />
        {props.children}
      </CardContent>
    </StyledCard>
  );
};

const StyledCard = styled(Card)`
  padding: 5px;
  margin-bottom: 10px;
`;

const Spacer = styled.div`
  margin-bottom: 10px;
`;
