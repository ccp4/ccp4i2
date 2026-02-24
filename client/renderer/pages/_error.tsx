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
