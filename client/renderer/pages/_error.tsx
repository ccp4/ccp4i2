/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
// Minimal Pages Router error page - overrides auto-generated _error.js
// which pulls in useContext-dependent code that fails during prerendering.
function ErrorPage() {
  return (
    <div style={{ padding: "2rem", textAlign: "center" }}>
      <h1>Error</h1>
      <p>An error occurred.</p>
    </div>
  );
}

ErrorPage.getInitialProps = () => {
  return { statusCode: 500 };
};

export default ErrorPage;
