import { serviceTypes } from "./experimentConstants";

export const handlerName = actionName => {
  actionName = actionName.replace(/_\d$/, "");
  return toUpperCase(removeUnderlines(serviceTypes[actionName]));
};

export const removeUnderlines = string => string.replace(/_/g, " ");

export const toUpperCase = string =>
  string.substr(0, 1).toUpperCase() + string.substr(1);
